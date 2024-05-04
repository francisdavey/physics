---
title: Welcome page
layout: default
---

This is it!

***
{% for post in site.posts %}
<h2 class="post-title"><a href="{{ post.permalink }}">{{ post.title }}</a> {{ post.date | date: "%e %B %Y" }}</h2>
{% if post.subtitle %}
<h3 class="post-subtitle">{{post.subtitle}}</h3>
{% endif %}
<p>{{ post.content }}</p>
{% endfor %}
